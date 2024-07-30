import os
import functools
import pathlib
import yaml
import gc
import pickle
import glob

import typing
import tqdm
import click
from click_option_group import optgroup
import yaml
import pytorch_lightning as pl
import torch
from pytorch_lightning.loggers import MLFlowLogger
# from pytorch_lightning.profiler import AdvancedProfiler
from pytorch_lightning.callbacks import TQDMProgressBar, EarlyStopping, ModelCheckpoint


@click.command()
@click.option(
    "--model-config-file",
    "model_config_files",
    help="The path to a YAML configuration file for the model.",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    default=tuple(),
    multiple=True,
)
@click.option(
    "--data-config-file",
    "data_config_files",
    help="The path to a YAML configuration file for the data.",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    default=tuple(),
    multiple=True,
)
@click.option(
    "--data-cache-directory",
    help="Path to cached data",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    required=True,
)
@click.option(
    "--output-directory",
    help="The path to an output directory",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    default=".",
    show_default=True,
)
@click.option(
    "--n-epochs",
    help="Number of epochs",
    type=int,
    default=200,
    show_default=True,
)
@click.option(
    "--n-gpus",
    help="Number of gpus",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "--learning-rate",
    help="Adam optimizer learning rate",
    type=float,
    default=0.001,
    show_default=True,
)
@click.option(
    "--n-processes",
    help="Number of processes",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--starting-model",
    help="Path to a starting model to load weights from",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    default=None,
)
def train_model(
    data_cache_directory: str,
    model_config_files: typing.Tuple[str, ...] = tuple(),
    data_config_files: typing.Tuple[str, ...] = tuple(),
    output_directory: str = ".",
    n_gpus: int = 1,
    n_epochs: int = 200,
    learning_rate: float = 0.001,
    n_processes: int = 4,
    starting_model: str = None
):
    from openff.nagl.config.training import TrainingConfig

    from openff.nagl.training.training import TrainingGNNModel, GNNModel, DGLMoleculeDataModule

    from pytorch_lightning.callbacks import ModelCheckpoint

    import torch
    torch.cuda.empty_cache()



    model_config_kwargs = {}
    for model_file in tqdm.tqdm(model_config_files, desc="Loading model configs"):
        with open(model_file, "r") as f:
            model_config_kwargs.update(yaml.load(f, Loader=yaml.FullLoader))
    
    data_config_kwargs = {}
    for data_file in tqdm.tqdm(data_config_files, desc="Loading data configs"):
        with open(data_file, "r") as f:
            data_config_kwargs.update(yaml.load(f, Loader=yaml.FullLoader))

    for v in data_config_kwargs.values():
        v["cache_directory"] = os.path.abspath(data_cache_directory)

    optimizer_config_kwargs = {
        "optimizer": "Adam",
        "learning_rate": learning_rate
    }

    training_config = TrainingConfig(
        model=model_config_kwargs,
        data=data_config_kwargs,
        optimizer=optimizer_config_kwargs,
    )

    training_model = TrainingGNNModel(training_config)
    # load weights from starting model
    if starting_model is not None:
        training_model.model = GNNModel.load(starting_model, eval_mode=False)

    # training_model.to("cuda")
    print(f"""
    Training model:
    {training_model}
    """)
    data_module = DGLMoleculeDataModule(
        training_config,
        n_processes=n_processes
    )
    data_module.prepare_data()
    data_module.setup()
    train_dataloader = data_module.train_dataloader()
    print(
        len(train_dataloader),
        len(train_dataloader.dataset),
        len(data_module._datasets["train"]),
        len(data_module._datasets["train"].datasets),
        data_module._dataset_configs["train"]
    )

    val_dataloader = data_module.val_dataloader()
    print(
        len(val_dataloader),
        len(val_dataloader.dataset),
        len(data_module._datasets["val"]),
        len(data_module._datasets["val"].datasets),
        data_module._dataset_configs["val"]
    )

    output_directory = pathlib.Path(output_directory)
    log_output_directory = output_directory / "mlruns"
    log_output_directory.mkdir(parents=True, exist_ok=True)
    logger = MLFlowLogger(
        experiment_name="defaults",
        save_dir=str(log_output_directory)
    )

    print(f"""
    Data:
    {data_module}

    Logger:
    {logger}
    -------------------------
    """
    )

    checkpoint_directory = output_directory / "checkpoints"
    last_checkpoint = checkpoint_directory / "last.ckpt"
    if not last_checkpoint.exists():
        last_checkpoint = None
    else:
        last_checkpoint = str(last_checkpoint)

    checkpoint_callback = ModelCheckpoint(
        save_last=True,
        dirpath=checkpoint_directory,
        monitor="val/loss",
    )
    callbacks = [
        TQDMProgressBar(),
        checkpoint_callback
    ]

    
    print(f"""
    Using {n_gpus} GPUs with {n_epochs} epochs
    Starting from {last_checkpoint}
    """)

    trainer = pl.Trainer(
        gpus=n_gpus,
        min_epochs=n_epochs,
        max_epochs=n_epochs,
        logger=logger,
        log_every_n_steps=1,
        check_val_every_n_epoch=1,
        enable_checkpointing=True,
        callbacks=callbacks,
        # profiler=profiler,
    )
    print(trainer)

    trainer.fit(
        training_model,
        datamodule=data_module,
        ckpt_path=last_checkpoint
    )

    if data_module.config.data.test.sources:
        trainer.test(
            training_model,
            datamodule=data_module,
        )

    metrics_file = str(output_directory / "metrics.pkl")
    with open(metrics_file, "wb") as f:
        metrics = (trainer.callback_metrics, trainer.logged_metrics)
        pickle.dump(metrics, f)
    print(f"Wrote metrics to {metrics_file}")

    print(f"{checkpoint_callback.best_model_path} : {checkpoint_callback.best_model_score}")

    best_model = TrainingGNNModel.load_from_checkpoint(
        str(checkpoint_callback.best_model_path)
    )

    best_training_model_path = str(output_directory / "best_training_model.pt")
    best_model_path = str(output_directory / "best_model.pt")

    torch.save(best_model, best_training_model_path)
    best_model.model.save(best_model_path)

    print(f"Saved best model to {best_model_path} (GNNModel.save)")
    print(f"Saved best training model to {best_training_model_path} (torch.save)")

    best_model_state_path = str(output_directory / "best_model_state_dict.pt")
    torch.save(best_model.state_dict(), best_model_state_path)
    print(f"Saved best training model state dict to {best_model_state_path}")




if __name__ == "__main__":
    train_model()
